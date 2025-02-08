# Открыть изображение как бинарный файл
with open("bad10.jpg", "rb") as file:
    data = file.read()

# Вывести первые 100 байт
print(data[:100])
